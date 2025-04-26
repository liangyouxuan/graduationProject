<template>
  <div class="home">
    <!-- 顶部导航栏 -->
    <el-header class="header">
      <div class="nav-container">
        <div class="logo">
          高分子相容性预测系统
        </div>
        <div class="auth-buttons">
          <template v-if="!isLoggedIn">
            <el-button type="success" class="custom-button" @click="$router.push('/login')">登录</el-button>
            <el-button type="success" class="custom-button" @click="$router.push('/register')">注册</el-button>
          </template>
          <template v-else>
            <el-dropdown>
              <el-button type="success" class="custom-button">
                {{ username }}<el-icon class="el-icon--right"><CaretBottom /></el-icon>
              </el-button>
              <template #dropdown>
                <el-dropdown-menu>
                  <el-dropdown-item @click="$router.push('/profile')">个人中心</el-dropdown-item>
                  <el-dropdown-item @click="handleLogout">退出登录</el-dropdown-item>
                </el-dropdown-menu>
              </template>
            </el-dropdown>
          </template>
        </div>
      </div>
    </el-header>

    <div class="main-container">
      <!-- 左侧导航栏 -->
      <div class="side-nav">
        <router-link to="/" class="nav-item" exact>
          <el-icon><House /></el-icon>
          <span>首页</span>
        </router-link>
        <router-link to="/prediction" class="nav-item">
          <el-icon><Search /></el-icon>
          <span>溶剂预测</span>
        </router-link>
      </div>

      <!-- 主要内容区域 -->
      <div class="content-area">
        <router-view></router-view>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import { useRouter } from 'vue-router'
import { ElMessage } from 'element-plus'
import { CaretBottom, House, Search } from '@element-plus/icons-vue'
import axios from 'axios'

const router = useRouter()
const isLoggedIn = ref(false)
const username = ref('')

// 获取用户信息
const getUserInfo = async () => {
  try {
    const token = localStorage.getItem('token')
    if (!token) return
    
    const response = await axios.get('http://localhost:8000/api/user/info/', {
      headers: {
        'Authorization': `Bearer ${token}`
      }
    })
    username.value = response.data.username
  } catch (error) {
    console.error('获取用户信息失败:', error)
  }
}

// 检查登录状态
const checkLoginStatus = async () => {
  const token = localStorage.getItem('token')
  isLoggedIn.value = !!token
  if (token) {
    try {
      await getUserInfo()
      axios.defaults.headers.common['Authorization'] = `Bearer ${token}`
    } catch (error) {
      if (error.response?.status === 401) {
        localStorage.removeItem('token')
        isLoggedIn.value = false
        router.push('/login')
      }
    }
  }
}

// 处理退出登录
const handleLogout = () => {
  localStorage.removeItem('token')
  delete axios.defaults.headers.common['Authorization']
  isLoggedIn.value = false
  ElMessage.success('已退出登录')
  router.push('/login')
}

onMounted(async () => {
  await checkLoginStatus()
})
</script>

<style scoped>
.home {
  min-height: 100vh;
  display: flex;
  flex-direction: column;
}

.header {
  background: #2c785c;
  padding: 0;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.nav-container {
  height: 60px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0 20px;
}

.logo {
  color: white;
  font-size: 24px;
  font-weight: bold;
  text-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);
  letter-spacing: 1px;
}

.auth-buttons {
  display: flex;
  gap: 10px;
}

.custom-button {
  background-color: #84fab0 !important;
  border-color: #84fab0 !important;
  color: #2c785c !important;
}

.custom-button:hover {
  background-color: white !important;
  border-color: white !important;
  color: #2c785c !important;
}

.main-container {
  flex: 1;
  display: flex;
  background-color: #f5f5f5;
}

.side-nav {
  width: 200px;
  background-color: white;
  padding: 20px 0;
  box-shadow: 2px 0 4px rgba(0, 0, 0, 0.1);
}

.nav-item {
  display: flex;
  align-items: center;
  padding: 12px 24px;
  color: #333;
  text-decoration: none;
  transition: all 0.3s ease;
}

.nav-item:hover {
  background-color: #f0f9eb;
  color: #2c785c;
}

.nav-item.router-link-active {
  background-color: #f0f9eb;
  color: #2c785c;
  border-right: 3px solid #2c785c;
}

.nav-item .el-icon {
  margin-right: 12px;
  font-size: 18px;
}

.content-area {
  flex: 1;
  padding: 20px;
  overflow-y: auto;
}
</style> 